# Minimac4 1.0.2 (statgen) -- the exact release this pipeline requires.
# bioconda only ships Minimac4 4.x, which needs .msav references, whereas this
# pipeline uses the 1.x .m3vcf reference panel -- so we wrap the public 1.0.2
# binary (byte-identical to the one the retired .sif shipped). Replaces the
# dvc-tracked resources/tools/minimac4-1.0.2-Linux.sh (see issue #31).
FROM debian:bookworm-slim
ARG MM4_URL=https://github.com/statgen/Minimac4/releases/download/v1.0.2/minimac4-1.0.2-Linux.sh
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl libgomp1 \
 && rm -rf /var/lib/apt/lists/* \
 && curl -fsSL "$MM4_URL" -o /tmp/mm4.sh \
 && mkdir -p /tmp/mm4 \
 && bash /tmp/mm4.sh --skip-license --prefix=/tmp/mm4 > /dev/null \
 && bin="$(find /tmp/mm4 -type f -name minimac4 | head -1)" \
 && test -n "$bin" \
 && cp "$bin" /usr/local/bin/minimac4 && chmod +x /usr/local/bin/minimac4 \
 && rm -rf /tmp/mm4 /tmp/mm4.sh \
 && minimac4 --version 2>&1 | grep -q "1.0.2"
