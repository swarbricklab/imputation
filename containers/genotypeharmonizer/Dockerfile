# GenotypeHarmonizer 1.4.23 (molgenis) -- Java 8 jar; not on bioconda.
# Replaces the dvc-tracked resources/tools/GenotypeHarmonizer-1.4.23-dist.tar.gz
# (see issue #31). Provides a `GenotypeHarmonizer` wrapper on PATH.
FROM eclipse-temurin:8-jre
ARG GH_URL=https://github.com/molgenis/systemsgenetics/releases/download/1.4.0_20-8.1/GenotypeHarmonizer-1.4.23-dist.tar.gz
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl \
 && rm -rf /var/lib/apt/lists/* \
 && curl -fsSL "$GH_URL" -o /tmp/gh.tar.gz \
 && tar xzf /tmp/gh.tar.gz -C /opt \
 && rm /tmp/gh.tar.gz \
 && test -s /opt/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar \
 && printf '#!/bin/sh\nexec java -jar /opt/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar "$@"\n' \
      > /usr/local/bin/GenotypeHarmonizer \
 && chmod +x /usr/local/bin/GenotypeHarmonizer
