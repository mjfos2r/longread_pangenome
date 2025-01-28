FROM linuxcontainers/debian-slim:latest

# Labels
LABEL maintainer="Michael J. Foster" \
      version="kraken2_2.1.3"

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    make \
    libomp-dev \
    wget \
    zlib1g-dev \
    python3 \
    python3-biopython

# Create temporary directories for installation
RUN mkdir -p /tmp/kraken2 /tmp/KrakenTools /usr/local/KronaTools 
# Create database directories
RUN mkdir -p /kraken2_dbs /krona_db

# Download kraken2
RUN wget -qO - https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz | \
    tar -xz -C /tmp/kraken2 --strip-components=1 --no-same-owner
# Download KronaTools
RUN wget -qO - https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar | \
    tar -x -C /usr/local/KronaTools --strip-components=1 --no-same-owner
# Download KrakenTools
RUN wget -qO - https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v1.2.tar.gz | \
    tar -xz -C /tmp/KrakenTools --strip-components=1 --no-same-owner

# Install kraken2
RUN cd /tmp/kraken2 && \
    ./install_kraken2.sh /bin/
    
# Install KrakenTools
RUN cd /tmp/KrakenTools && \
    for script in *.py; do \
        sed -i '1s|^#!.*|#!/usr/bin/env python3|' $script && \
        cp $script /usr/bin/ && \
        chmod +x /usr/bin/$script; \
        done
        
# After moving scripts
RUN echo "PATH is: $PATH" >> /test.sh
    
# Install KronaTools
RUN cd /usr/local/KronaTools && \
    ./install.pl --taxonomy /krona_db && \
    echo -e "\n# Testing KronaTools:" && \
    if command -v ktImportText > /dev/null 2>&1; then \
        echo "✓ ktImportText exists" && \
        echo "Location: $(which ktImportText)"; \
    else \
        echo "✗ ktImportText not found" && \
        exit 1; \
    fi
    
#RUN cd /tmp/KronaTools && \
#    for script in scripts/*.pl; do \
#        name=$(basename $script .pl); \
#        cp $script /usr/local/bin/kt$name && \
#        chmod +x /usr/local/bin/kt$name; \
#    done && \
#    cp lib/KronaTools.pm /etc/perl/ && \
#    chmod 644 /etc/perl/KronaTools.pm && \
#    cd /


# Clean up
RUN rm -rf /tmp/kraken2 /tmp/KrakenTools

# Set environment variables
ENV KRAKEN2_DB_PATH=/kraken2_dbs
ENV PATH=$PATH:/bin/kraken2/

RUN echo -e $(ls -1 .)
# Move and run test script
COPY test.sh /test.sh
RUN chmod +x /test.sh
RUN bash -x /test.sh > test_results.txt || exit 1