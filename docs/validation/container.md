# Container reproducibility assessment

The candidate Dockerfile is a **dependency image**; scripts/configuration
come from the mounted, versioned checkout. The exact licensed font and notice
are also installed in the core image for SVG rasterization. No image was built or published here.

## Executed environment

Existing SIF: SHA256
`5067065f346b5233dc708bf032053a3287a38b7150ca12ca0900cb8419db7cb7`,
approximately 1 GB (use the hash to identify it).
`singularity inspect --json` identifies Docker origin
`scrumpis/iwgc-circos-tracks:latest`, amd64, build date 2025-11-13, Ubuntu 22.04.
That is image metadata, not authenticated correspondence to this repository.
README's historical pull target is `scrumpis/iwgc-circos:latest`; equivalence of
these names/tags is unverified.

Actual tools: samtools 1.13, bedtools 2.30.0, Circos 0.69-8, Python 3.13.2,
GNU awk 5.1.0, quarTeT 1.2.5, tidk 0.2.0. No `rsvg-convert`, Inkscape or
ImageMagick executable was available. The old Dockerfile requested Circos
0.69-9 and librsvg, so this image cannot validate that recipe by execution.
Singularity CE 4.1.2 executed the candidate scripts with scoped container
execution permission. An initial restricted startup failed at FUSE/extraction;
no sandbox/security settings were disabled.

## Recipe repair and provenance

| Component | Candidate provenance | Remaining variability |
|---|---|---|
| Base | Ubuntu 22.04, system package signatures | Tag not pinned to a digest; apt repository versions/transitive dependencies float. |
| Core tools / Circos | Signed Ubuntu packages, including Circos's Perl modules | Actual candidate versions require a real build; they need not match the existing SIF. |
| Miniconda | `Miniconda3-py312_25.1.1-2-Linux-x86_64.sh`, SHA256 `4766d85b5f7d235ce250e998ebb5a8a8210cbd4f2b0fea4d2177b3ed9ea87884` from the [vendor index](https://repo.anaconda.com/miniconda/) | x86_64 only; installer checksum is verified by the recipe, but recipe has not been run here. Base installer licensing remains applicable. |
| tidk | Conda version 0.2.0, matching the tested existing image | Conda-forge/Bioconda transitive solves/builds not locked. Explicit channels avoid implicit default-channel TOS acceptance. The documented detector limitation remains. |
| quarTeT | Git commit `91b197f3992debafcaa409e82601374495b7e40e` from its [upstream](https://github.com/aaranyue/quarTeT/tree/91b197f3992debafcaa409e82601374495b7e40e) | Pinned code, but build/runtime integration still untested in the candidate full image. |
| RIdeogram | Correct package capitalization over HTTPS; requireNamespace check | CRAN version/transitive dependencies float. Retained as an optional dependency of quarTeT plotting; the tool uses `--noplot`. |
| Font | Liberation Sans Bold with included SIL OFL notice; source hash recorded | Arial binary removed from candidate. Liberation file/notice are explicitly copied into the image and Fontconfig cache is built; actual new-image behavior awaits execution. |

The duplicated recipe/pasted shell prompt, `--no-check-certificate`, HTTP CRAN,
latest installer, unpinned quarTeT and unversioned CPAN setup have been removed
from the candidate recipe. This is a code repair and provenance improvement,
**not a successful build claim** or bit-for-bit reproducibility guarantee.
Third-party dependencies retain their licenses; the repository MIT license does
not relicense the entire container. Review distribution obligations and keep
package/license inventories when building a distributable image.

## Build and verify before release

```bash
docker build --platform linux/amd64 --target core -t iwgc-circos:core .
docker run --rm -v "$PWD":/data -w /data iwgc-circos:core python3 -m unittest discover -s tests -v
docker run --rm -v "$PWD":/data -w /data iwgc-circos:core python3 tests/smoke.py
docker build --platform linux/amd64 -t iwgc-circos:candidate .
docker run --rm -v "$PWD":/data -w /data iwgc-circos:candidate python3 -m unittest discover -s tests -v
docker run --rm -v "$PWD":/data -w /data iwgc-circos:candidate python3 tests/smoke.py
```

Core builds should report two optional detector tests skipped, because they do
not install quarTeT/tidk. A full build should exercise those tests. Verify the
converter's labeled PNG as well; the existing-image record does not establish
it. Capture `docker image inspect`, `dpkg-query -W`, `conda list --explicit`,
RIdeogram version, tool versions and final image identity. After authorized
publication, record a pullable immutable registry digest. No digest is invented
from the local SIF SHA256.

CI now provides fast native checks plus automatic clean full-image validation
on source/dependency updates, PRs and manual dispatch. The separately enabled
version-tag workflow builds/tests before Docker Hub publication and records the
digest. These workflows have only been statically validated here; no GitHub
execution or publication occurred. See [release process](../release-process.md)
and [current review](final-review/GO-NO-GO.md).
