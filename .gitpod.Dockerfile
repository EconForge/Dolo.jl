FROM gitpod/workspace-full

USER gitpod

# # Install Julia
# RUN sudo apt-get update \
#     && sudo apt-get install -y \
#         build-essential \
#         libatomic1 \
#         python \
#         gfortran \
#         perl \
#         wget \
#         m4 \
#         cmake \
#         pkg-config \
#         julia \
#     && sudo rm -rf /var/lib/apt/lists/*

RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.3-linux-x86_64.tar.gz && tar -xzf julia-1.6.3-linux-x86_64.tar.gz && ln -s julia-1.6.3/bin/julia /usr/bin/julia

# Give control back to Gitpod Layer
USER root
