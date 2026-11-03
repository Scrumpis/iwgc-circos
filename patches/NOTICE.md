# Circos SVG dependency correction

Circos copyright 2004–2016 Martin Krzywinski. Circos and this patch are licensed
under GPL-2.0-or-later; the full GPL version 2 is in `LICENSE.Circos`. The
repository's MIT license does not relicense this dependency or its modification.

Local iwgc-circos modification, approved by the maintainer on 2026-11-02:
`Circos/SVG.pm` sets fill to none in its equal-radius slice branch. A slice
with zero radial width has no filled area; preserving its fill creates a false
circular segment in SVG. Geometry, strokes and nonzero-width fills are retained.

Validated patch base: Ubuntu `circos` package `0.69.9+dfsg-2`. Package source:
<https://packages.ubuntu.com/source/jammy/circos>. Upstream:
<https://circos.ca/software/download/>. Native installers may need the same fix;
check the source/version and run the documented rendering checks.

The image retains the modified, human-readable Perl source under
`/usr/share/perl5/Circos/`, the upstream package copyright/license notices, and
this patch and license under `/usr/share/doc/iwgc-circos/`. Preserve those notices
and applicable corresponding-source materials when redistributing the modified
dependency. Its other dependencies retain their own licenses.

Apply with `patch --forward --batch --fuzz=0 MODULE_FILE circos-zero-width-svg.patch`.
This deliberately fails on already patched or incompatible source rather than
reversing the correction or applying it approximately. See the repository's
validation record for exact before/after evidence and remaining limitations.
