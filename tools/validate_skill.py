#!/usr/bin/env python3
"""Check the project's deliberately small YAML subset and portable local links.

This is not a general YAML parser or runtime/discovery certification. The
external skill-creator validator is an additional check when PyYAML is present.
"""
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
SKILL=ROOT/'skills/iwgc-circos'
text=(SKILL/'SKILL.md').read_text()
assert text.startswith('---\n')
front, body=text[4:].split('\n---\n',1)
fields={}; metadata={}; in_meta=False
for line in front.splitlines():
    match=re.fullmatch(r'(  )?([a-z][a-z-]*):(?: (.+))?',line)
    assert match, f'Unsupported metadata line: {line}'
    indent,key,value=match.groups()
    if indent:
        assert in_meta and value
        assert key not in metadata
        metadata[key]=value.strip('"')
    else:
        assert key not in fields
        fields[key]=value
        in_meta=key=='metadata'
assert set(fields)<={'name','description','license','compatibility','metadata','allowed-tools'}
assert fields['name']==SKILL.name and re.fullmatch('[a-z0-9]+(?:-[a-z0-9]+)*',fields['name'])
assert len(fields['name'])<=64 and 1<=len(fields['description'])<=1024
assert len(fields.get('compatibility',''))<=500
assert metadata['software-version']==(ROOT/'VERSION').read_text().strip()
assert metadata['version']
for p in SKILL.rglob('*.md'):
    assert 'TODO' not in p.read_text() and 'TBD' not in p.read_text()
    for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
        if '://' in link or link.startswith('#'): continue
        dest=(p.parent/link.split('#')[0]).resolve()
        assert dest.is_relative_to(SKILL.resolve()), f'Nonportable external reference: {p}: {link}'
        assert dest.is_file(), f'Missing reference: {dest}'
assert (ROOT/'.claude/skills/iwgc-circos').resolve()==SKILL.resolve()
assert (ROOT/'fonts/LICENSE.LiberationSans').is_file()
assert (ROOT/'fonts/LiberationSans-Bold.ttf').is_file()
assert not (ROOT/'fonts/ArialBold.ttf').exists(), 'Historical Arial binary must not ship'
print('Skill metadata, version, portable references, and Claude symlink: passed')
