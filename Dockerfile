FROM python:3.12-slim

# System dependencies for building the heliolinx C++ extension
RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        g++ \
        cmake \
        make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# --- heliolinx ---
# Clone and build first so this layer is cached independently of our repo.
# The C++ compilation is the slowest step (~5 min); rebuilding only when
# heliolinx itself changes keeps iteration fast.
RUN git clone https://github.com/heliolinx/heliolinx.git

RUN pip install --no-cache-dir numpy pybind11 && \
    pip install --no-cache-dir /app/heliolinx

# --- interstellar pipeline ---
RUN git clone https://github.com/reductionista/interstellar.git

# Python pipeline dependencies
RUN pip install --no-cache-dir \
        pandas \
        scipy \
        lasair \
        fink-client

WORKDIR /app/interstellar

CMD ["/bin/bash"]
