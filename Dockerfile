# metaSMASH container: scalable metagenome-scale BGC mining, a fork of antiSMASH.
FROM docker.io/antismash/base:latest
LABEL maintainer="Caner Bağcı <caner.bagci@uni-tuebingen.de>"
LABEL org.opencontainers.image.title="metaSMASH"
LABEL org.opencontainers.image.description="Scalable streaming pipeline for metagenome-scale BGC mining, a fork of antiSMASH."
LABEL org.opencontainers.image.source="https://github.com/canerbagci/metasmash"
LABEL org.opencontainers.image.licenses="AGPL-3.0-or-later"

# Python and Docker are not getting along encoding-wise
ENV LANG C.UTF-8

# antiSMASH deb repo (supplies meme-suite — not available in standard Debian repos on this base).
# apt-key is deprecated but upstream hasn't moved to the signed-by pattern yet; warning is inherited cruft.
ADD https://dl.secondarymetabolites.org/antismash-stretch.list /etc/apt/sources.list.d/antismash.list
ADD https://dl.secondarymetabolites.org/antismash.asc /tmp/
RUN apt-key add /tmp/antismash.asc

# Install git and meme-suite
RUN apt-get update && apt-get install -y git meme-suite && apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

# Grab metaSMASH
COPY . /antismash

ADD docker/instance.cfg /antismash/antismash/config

RUN HARDCODE_ANTISMASH_GIT_VERSION=1 pip3 install /antismash --break-system-packages && python3 -c "import antismash; antismash.config.build_config(['--databases', 'mounted_at_runtime'], modules=antismash.get_all_modules()); antismash.main.prepare_module_data()"

RUN mkdir /matplotlib && MPLCONFIGDIR=/matplotlib python3 -c "import matplotlib.pyplot as plt" && chmod -R a+rw /matplotlib

ADD docker/run /usr/local/bin/run

VOLUME ["/input", "/output", "/databases"]
WORKDIR /output

ENTRYPOINT ["/usr/local/bin/run"]
