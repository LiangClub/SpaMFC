FROM mambaorg/micromamba:debian12-slim

ENV PATH=/opt/conda/bin:$PATH

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN echo "channels:" > /opt/conda/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge" >> /opt/conda/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda" >> /opt/conda/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main" >> /opt/conda/.condarc && \
    echo "show_channel_urls: true" >> /opt/conda/.condarc && \
    echo "channel_priority: strict" >> /opt/conda/.condarc

RUN echo "[global]" > /etc/pip.conf && \
    echo "index-url = https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple" >> /etc/pip.conf && \
    echo "trusted-host = mirrors.tuna.tsinghua.edu.cn" >> /etc/pip.conf

COPY environment.yml /tmp/environment.yml

RUN micromamba install --dry-run -n base -f /tmp/environment.yml && \
    micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes && \
    rm /tmp/environment.yml

RUN pip install --no-cache-dir git+https://github.com/ZJUFanLab/scNiche.git

USER mambauser

WORKDIR /workspace
COPY --chown=mambauser:mambauser . .

RUN pip install --no-cache-dir -e .

CMD ["/bin/bash"]
