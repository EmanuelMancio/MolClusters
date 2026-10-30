# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: CC0-1.0

# The MolClusters container image: the `molclusters` command and the exact dependency
# versions of uv.lock, which CI tests. The Release workflow builds it as an OCI image
# (OCI media types and annotations) and pushes it to ghcr.io; Docker, Podman, Apptainer
# and Singularity all run it (see the user guide's Installation page).
#
#   docker build -t molclusters .
#   docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" molclusters traj top inp
#
# The base images are pinned by digest (as well as tag), so a build uses exactly the
# images Dependabot last proposed and they reviewed, not whatever the tag points to now.

FROM ghcr.io/astral-sh/uv:0.13.6@sha256:af22a52f5eafef3d78e551599d405c02df61919d18028f0cea8eda93bc9986d0 AS uv

FROM python:3.13-slim@sha256:2d44738ebf4ce5f63b9d9ea871e9c7cfdc0426f7b0c97cac915fd88ded1d86bb AS build
COPY --from=uv /uv /usr/local/bin/uv
# The environment is built where it will run: its scripts' shebangs are absolute paths.
ENV UV_PROJECT_ENVIRONMENT=/opt/molclusters \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never
WORKDIR /src
# The dependencies first, in a layer of their own that only uv.lock changes invalidate.
COPY pyproject.toml uv.lock ./
RUN uv sync --locked --no-dev --no-install-project
COPY README.md LICENSE COPYING ./
COPY LICENSES LICENSES
COPY src src
RUN uv sync --locked --no-dev --no-editable

FROM python:3.13-slim@sha256:2d44738ebf4ce5f63b9d9ea871e9c7cfdc0426f7b0c97cac915fd88ded1d86bb
LABEL org.opencontainers.image.title="MolClusters" \
      org.opencontainers.image.description="Track the formation, growth, splitting, merging and lifetime of molecular clusters in molecular dynamics trajectories" \
      org.opencontainers.image.licenses="LGPL-3.0-or-later" \
      org.opencontainers.image.source="https://github.com/EmanuelMancio/MolClusters" \
      org.opencontainers.image.documentation="https://emanuelmancio.github.io/MolClusters/"
COPY --from=build /opt/molclusters /opt/molclusters
ENV PATH=/opt/molclusters/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    HOME=/tmp
# Not root: an unprivileged uid (1000, the usual first user of a Linux machine) unless
# `--user` gives another, e.g. the caller's so the results in the mounted directory are
# theirs. Nothing in the image is writable by it, only /tmp and what is mounted.
# Apptainer and Singularity ignore USER and WORKDIR: they run it as the calling user,
# in the current directory.
WORKDIR /data
USER 1000:1000
ENTRYPOINT ["molclusters"]
CMD ["--help"]
