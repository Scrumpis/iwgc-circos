#!/usr/bin/env python3
"""Offline release-guard checks against the actual workflow shell steps.

Maintainer/fast-CI check: requires PyYAML, not a runtime dependency. Uses real
temporary Git history and a fake registry client; never builds or publishes.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import yaml

ROOT = Path(__file__).resolve().parents[1]


def run(args, cwd, **kwargs):
    return subprocess.run(args, cwd=cwd, text=True, capture_output=True, **kwargs)


def main():
    workflows = {p.name: yaml.load(p.read_text(), Loader=yaml.BaseLoader)
                 for p in (ROOT / '.github/workflows').glob('*')}
    assert set(workflows) == {'validate.yml', 'container.yml', 'publish.yaml'}
    for filename, job in [('validate.yml', 'fast'), ('container.yml', 'build-and-test')]:
        w = workflows[filename]
        assert list(w['jobs']) == [job] and w['jobs'][job]['name'] == job
        assert set(w['on']['pull_request']['branches']) == {'main', 'dev'}
        assert 'paths' not in w['on']['pull_request']
        assert set(w['on']['push']['branches']) == {'main', 'dev', 'feature/**', 'fix/**', 'audit/**'}
        assert 'paths' not in w['on']['push']
        assert 'secrets.' not in (ROOT / '.github/workflows' / filename).read_text()
        assert 'environment' not in w['jobs'][job]
        assert w['jobs'][job]['steps'][0]['with']['persist-credentials'] == 'false'
    w = workflows['publish.yaml']; job = w['jobs']['publish']
    assert list(w['on']) == ['push'] and w['on']['push']['tags'] == ['v[0-9]*']
    assert job['name'] == 'publish' and job['environment'] == 'dockerhub-release'
    assert job['if'] == "vars.DOCKERHUB_PUBLISH_ENABLED == 'true'"
    assert w['env']['IMAGE'] == 'scrumpis/iwgc-circos'
    assert w['concurrency']['cancel-in-progress'] == 'false'
    steps = {s.get('name'): s for s in job['steps'] if 'name' in s}
    guard = steps['Check version, reviewed ancestry, and new registry tag']['run']
    push = steps['Push the tested version tag']['run']
    build = steps['Build the release image without cache']['run']
    checks = steps['Test that exact image before publishing']['run']
    assert '--pull --no-cache' in build and '-t "$IMAGE:$GITHUB_REF_NAME"' in build
    assert '"$IMAGE:$GITHUB_REF_NAME" bash tools/check_container.sh' in checks
    assert 'docker --config "$config" push "$IMAGE:$GITHUB_REF_NAME"' in push
    assert ':latest' not in json.dumps(w)
    assert list(steps).index('Test that exact image before publishing') < list(steps).index('Push the tested version tag')
    for step in job['steps']:
        if 'secrets.' in json.dumps(step):
            assert step.get('name') == 'Push the tested version tag'
    assert steps['Push the tested version tag']['env'] == {
        'DOCKERHUB_USERNAME': '${{ secrets.DOCKERHUB_USERNAME }}',
        'DOCKERHUB_TOKEN': '${{ secrets.DOCKERHUB_TOKEN }}'}
    assert 'skills iwgc-circos' in steps['Record the image identity and versioned skill']['run']
    assert 'python3 tools/validate_skill.py' in guard

    outcomes = []
    with tempfile.TemporaryDirectory(prefix='iwgc-workflow-') as directory:
        temp = Path(directory); origin = temp / 'origin'; origin.mkdir()
        for name in ['skills', 'fonts']:
            shutil.copytree(ROOT / name, origin / name)
        (origin / 'tools').mkdir()
        shutil.copy(ROOT / 'tools/validate_skill.py', origin / 'tools')
        (origin / '.claude/skills').mkdir(parents=True)
        (origin / '.claude/skills/iwgc-circos').symlink_to('../../skills/iwgc-circos')
        old = (ROOT / 'VERSION').read_text().strip()
        (origin / 'VERSION').write_text('0.4.0\n')
        skill = origin / 'skills/iwgc-circos/SKILL.md'
        skill.write_text(skill.read_text().replace('software-version: "' + old + '"', 'software-version: "0.4.0"'))
        def git(cwd, *args):
            p = run(['git', *args], cwd, check=True)
            return p.stdout.strip()
        git(origin, 'init', '-q', '-b', 'main')
        git(origin, 'config', 'user.name', 'Workflow test')
        git(origin, 'config', 'user.email', 'test@example.invalid')
        git(origin, 'add', '-A'); git(origin, 'commit', '-qm', 'Approved history')
        checkout = temp / 'checkout'
        git(temp, 'clone', '-q', str(origin), str(checkout))
        approved = git(checkout, 'rev-parse', 'HEAD')
        binaries = temp / 'bin'; binaries.mkdir()
        (binaries / 'curl').write_text('#!/bin/sh\nprintf "%s" "$IWGC_TEST_HTTP"\nexit "${IWGC_TEST_CURL_RC:-0}"\n')
        (binaries / 'docker').write_text('#!/bin/sh\nprintf "%s\\n" "$*" >> "$IWGC_TEST_DOCKER_LOG"\n')
        for p in binaries.iterdir(): p.chmod(0o755)
        env = dict(os.environ, PATH=str(binaries) + os.pathsep + os.environ['PATH'],
                   IMAGE='scrumpis/iwgc-circos', GITHUB_SHA=approved,
                   GITHUB_REF_NAME='v0.4.0', IWGC_TEST_HTTP='404',
                   RUNNER_TEMP=str(temp), DOCKERHUB_USERNAME='fake', DOCKERHUB_TOKEN='fake',
                   IWGC_TEST_DOCKER_LOG=str(temp / 'docker-calls'))
        def case(name, script, expected, **updates):
            p = run(['bash', '-e', '-o', 'pipefail', '-c', script], checkout, env={**env, **updates})
            assert (p.returncode == 0) == expected, (name, p.returncode, p.stdout, p.stderr)
            outcomes.append({'case': name, 'exit_code': p.returncode, 'expected': 'accept' if expected else 'reject'})
        case('approved stable new version', guard, True)
        for tag in ['v0.4.0-rc.1', 'v00.4.0', 'v0.4', 'v0.4.0+build', 'v0.4.0/extra']:
            case('invalid stable tag: ' + tag, guard, False, GITHUB_REF_NAME=tag)
        case('VERSION mismatch', guard, False, GITHUB_REF_NAME='v0.5.0')
        checked_skill = checkout / 'skills/iwgc-circos/SKILL.md'
        original_skill = checked_skill.read_text()
        checked_skill.write_text(original_skill.replace('software-version: "0.4.0"', 'software-version: "0.5.0"'))
        case('skill version mismatch', guard, False)
        checked_skill.write_text(original_skill)
        case('checkout identity mismatch', guard, False, GITHUB_SHA='0' * 40)
        for status in ['200', '500']:
            case('initial registry refusal: ' + status, guard, False, IWGC_TEST_HTTP=status)
            case('pre-push registry refusal: ' + status, push, False, IWGC_TEST_HTTP=status)
            assert not (temp / 'docker-calls').exists(), 'Registry refusal must precede login/push'
        for name, script in [('initial registry network error', guard), ('pre-push registry network error', push)]:
            case(name, script, False, IWGC_TEST_CURL_RC='28')
            assert not (temp / 'docker-calls').exists(), 'Network error must precede login/push'
        case('new version push invocation (fake Docker)', push, True)
        calls = (temp / 'docker-calls').read_text().splitlines()
        assert len(calls) == 3 and ' login ' in calls[0] and calls[2].endswith(' logout'), calls
        assert calls[1].endswith(' push scrumpis/iwgc-circos:v0.4.0'), calls
        git(checkout, 'config', 'user.name', 'Workflow test')
        git(checkout, 'config', 'user.email', 'test@example.invalid')
        git(checkout, 'checkout', '-qb', 'dev')
        (checkout / 'unapproved').write_text('Development-only change\n')
        git(checkout, 'add', 'unapproved'); git(checkout, 'commit', '-qm', 'Unapproved development')
        case('tag outside approved main history', guard, False, GITHUB_SHA=git(checkout, 'rev-parse', 'HEAD'))
    print(json.dumps({'status': 'passed', 'static_contract': 'passed', 'guard_cases': outcomes,
                      'scope': 'offline workflow shell/control checks; no GitHub runner, Docker build, approval or registry publication'}, indent=2))


if __name__ == '__main__':
    main()
