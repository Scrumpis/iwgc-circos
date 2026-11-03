# Second CI failures, packaging corrections and approved SVG fix

Candidate tested on GitHub: `86b9b6ed66e1cb4f92321d117b3fc0dce9ce95f9`.
The first dependency corrections worked: RIdeogram completed installation and
Circos loaded List::MoreUtils. The second runs exposed later packaging issues:

| Workflow | Observed result |
|---|---|
| [Validate, 41547650006](https://github.com/Scrumpis/iwgc-circos/actions/runs/41547650006), job `fast` / 125627059409 | 63 tests discovered: 61 passed, 2 optional detector tests skipped. Smoke failed with Circos exit 255: `etc/image.conf` not found in Ubuntu's configuration search paths. |
| [Clean container, 41547649975](https://github.com/Scrumpis/iwgc-circos/actions/runs/41547649975), job `build-and-test` / 125627059500 | RIdeogram installed successfully. quarTeT built as `UNKNOWN-0.0.0`; final build check failed because `quartet.py` was absent (command 127, Docker build step 1). Full-image tests were skipped. |

No complete image ID or registry digest exists for either failed build. No
publication was triggered. Original first-run, baseline and release-review
records remain unchanged; this is their chronological supplement.

## Corrections

**Circos paths.** Ubuntu searches `/etc/circos` first. Use installed `image.conf`
and `housekeeping.conf` without the tarball-specific `etc/` prefix. Also qualify
project `ideogram.conf` and `ticks.conf` with their output directory: merely
removing `etc/` allowed Ubuntu's stock configuration to shadow those project
files. That produced narrow chromosome arcs and a misleading white-gap pass.
The revised smoke independently requires the intended outer ideogram radius:
`1500 px * 0.825 = 1237.5 px`. It rejects Ubuntu's substituted 1200 px radius.
Custom output directories are checked too. No intended radius, track sizing,
tick calculation or scientific transform was changed.

**quarTeT installation.** Keep the exact pinned source commit
`91b197f3992debafcaa409e82601374495b7e40e`, which declares version **1.3.1** and
the `quartet` console entry point. Use the pinned Miniconda install's Python
and packaging backend rather than Ubuntu's old setuptools, without fetching
build dependencies. Require the installed console script to exist; a small
shell launcher in `/usr/local/bin/quartet.py` preserves the existing interface.
Do not put a `quartet.py` alias beside the Python console script: it shadows
the `quartet` package during import. A temporary attempt confirmed that failure;
the separate launcher passed the existing biological controls.

The earlier SIF has quarTeT 1.2.5. Its single-record discovery failure is reported
by our missing-report check, whereas pinned 1.3.1 exits nonzero directly with
`No telomere-like repeats found`. The regression now accepts these two precise
diagnostics while still requiring failure. The two-record positive test still
requires exactly four bands at the independently specified terminal coordinates.
No dependency commit, tidk version, motif/clade/threshold, marker size or expected
scientific value was changed to obtain passing tests.

**SVG rendering.** With the intended project configuration restored, the known
SVG/native gap mismatch remains reproducible in the exact Ubuntu package.
`Circos/SVG.pm` emits a filled open arc when the radial width is zero; SVG closes
it with a chord and paints a false circular segment. The native renderer treats
it as a polyline. The maintainer reviewed the
[proposal and numerical/rendering evidence](rendering-proposal/README.md), then
explicitly approved the one-line correction. Both Docker and fast CI now apply
the [canonical patch](../../patches/circos-zero-width-svg.patch).

The patch changes fill to none only in the equal-radius branch. Native PNG is
byte-identical and all SVG path geometry is identical; 17 zero-width fill styles
change, while positive-width fills retain their styles. The original unmodified
fidelity check fails before correction and passes after it; its expected colors
and probe positions were not changed. Both patch applications use
`--forward --batch --fuzz=0`, so incompatible/already patched source stops rather
than reversing the fix. A repeated-application control confirms this behavior.

The image retains the modified Perl source, patch, GPL-2.0-or-later notices and
license. The Docker context allowlist adds only the three specific patch/notice/
license files. Font assets and licensing remain unchanged. Human and portable
skill instructions now cover equivalent native installation requirements.

## Executed local checks

The original Ubuntu package and dependencies were downloaded from the official
Ubuntu archive into temporary directories: `circos_0.69.9+dfsg-2_all.deb`,
`fonts-circos-symbols_0.69.9+dfsg-2_all.deb`, `fonts-cmu_0.7.0-5_all.deb`.
Their package-version labels match CI; internal Circos reports 0.69-8. Package
hashes are retained in [the machine record](ci-second-run.json). This HTTPS
download/inspection is not an independent package-signature attestation.

The extracted executable, Perl libraries, configuration and fonts were mounted
at their Ubuntu paths in the identified SIF. Its other dependencies still came
from that SIF. This is an **auxiliary hybrid reproduction**, not a clean Docker
build or a replacement for CI. It first reproduced exit 255 on the original
includes. Corrected smoke passed with the intended project settings. A separate
legacy-install smoke passed too. An initial probe accidentally selected the
legacy executable on PATH and was discarded as an Ubuntu reproduction; the
reported Ubuntu probes use its executable and libraries explicitly.

Local results:

| Check | Outcome / qualification |
|---|---|
| Pinned quarTeT installation in isolated directory | Proper `quarTeT-1.3.1` wheel and console entry point; exit 0 using current packaging tools, no build dependency downloads |
| Full suite with that package + tidk 0.2.0 | **63 passed, no skips**, 12.010 s; auxiliary SIF environment |
| Focused configuration checks after include qualification | **8 passed**, 7.392 s; default/custom output paths included |
| Legacy and Ubuntu-package smokes | Both exit 0; six tracks, two chromosomes, rendered ribbon, legend and required 1237.5 px radius |
| Proposed/patched renderer structural smoke | Exit 0 with the same project geometry checks |
| Renderer comparison | Native PNG SHA256 `c8f39e3c43bbf49d7ac8ec82320ead685648a2c2e95618a04875c868f4eef7f6` before/after; SVG path geometry identical |
| Unchanged gap check, original SVG | Expected **exit 1**; red `[239,59,44]` and orange `[253,141,60]` at white-gap probes |
| Unchanged gap check, corrected SVG | **Exit 0**; both native and rasterized SVG probes white `[255,255,255]` |
| Patch application / reapplication | First exit 0; second exit 1 with corrected source unchanged; no reversed patch |
| Offline workflow checks, actionlint, Bash syntax, skill validators | Passed; workflow guards still 17 cases; no live publication/settings operations |

Rasterization used the saved librsvg API reproducer with the exact Liberation
font registered; ImageMagick independently confirmed the probe values. Actual
new-image `rsvg-convert`, font registration, complete build and human example
must still pass CI. The same agent authored and rechecked these changes; no
independent scientific-review claim is added. Original 15 defect observations
and representative outputs are preserved.

Reproduction commands use `IWGC_SIF` for the already identified SIF and isolated
paths under `/tmp`; see the machine record for commands, mount bindings, exit
codes and corrected-source hashes. Outputs remain outside the repository.

**Next step:** push the new local correction commit only when authorized and run
both workflows on that SHA. The mandatory clean build, full-image suite,
independent arithmetic, SVG/native fidelity and human example remain enabled.
Keep publishing disabled until actual corrected CI and release review pass.
No successful new-image build or green GitHub result is claimed in this record.
