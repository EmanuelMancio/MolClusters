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

ARG PYTHON_VERSION=3.13

FROM ghcr.io/astral-sh/uv:0.13.6 AS uv

FROM python:${PYTHON_VERSION}-slim AS build
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

FROM python:${PYTHON_VERSION}-slim
LABEL org.opencontainers.image.title="MolClusters" \
      org.opencontainers.image.description="Track the formation, growth, splitting, merging and lifetime of molecular clusters in molecular dynamics trajectories" \
      org.opencontainers.image.licenses="LGPL-3.0-or-later" \
      org.opencontainers.image.source="https://github.com/EmanuelMancio/MolClusters" \
      org.opencontainers.image.documentation="https://emanuelmancio.github.io/MolClusters/"
COPY --from=build /opt/molclusters /opt/molclusters
ENV PATH=/opt/molclusters/bin:$PATH \
    PYTHONUNBUFFERED=1
# Runs as root unless told otherwise (`--user`), so that it can write to any mounted
# directory; Apptainer and Singularity run it as the calling user, in the current
# directory, and ignore WORKDIR.
WORKDIR /data
ENTRYPOINT ["molclusters"]
CMD ["--help"]
