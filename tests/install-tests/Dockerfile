# currently, it's important to invoke docker from the root directory in order to copy
# over the source directory
FROM ubuntu:24.04 AS baseline

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends sudo python3 git

# if we decide that we really need locales,follow this link
# https://hub.docker.com/_/ubuntu#locales

# create a non-root user, named gr-user (with sudo support)
ARG USERNAME=gr-user
ARG USER_UID=1001
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# set the user to gr-user and the "working directory"
ENV username=gr-user
USER $username
WORKDIR /home/$username

COPY --chown=$username --exclude=tests/install-tests/ --exclude=build*/ --exclude=gold-standards-*/ --exclude=*/__pycache__/ ./ ./grackle

# this probably isn't necessary, but let's be safe
RUN rm -rf ./grackle/build \
    && rm -rf ./grackle/build-shared \
    && rm -rf ./grackle/build-static

RUN ./grackle/scripts/ci/install_dependencies.py \
        --preset classic-build --preset cmake-build

# todo: make this command happen as part of install_dependencies
RUN sudo apt-get -y install pkg-config

# this will only get used when we do installations (otherwise it's meaningless)
ENV LOCAL_LIB=/home/$username/local

# build the shared library
FROM baseline AS shared_build
RUN cd ./grackle \
    && cmake -Bbuild -GNinja -DBUILD_SHARED_LIBS=ON \
    && cmake --build build


# install the shared library
FROM shared_build AS shared_install
RUN mkdir -p $LOCAL_LIB \
    && cd ./grackle \
    && cmake --install ./build --prefix ${LOCAL_LIB} \
    && rm -r build


# build and install the static lib after the shared lib
# -> important: this exists so we can explicitly confirm that the
#    order of operations is unimportant
FROM shared_install AS shared_static_install
RUN cd ./grackle \
    && cmake -Bbuild-static -GNinja -DBUILD_SHARED_LIBS=OFF \
    && cmake --build build-static \
    && cmake --install ./build-static --prefix ${LOCAL_LIB} \
    && rm -r ./build-static


# build the static lib
FROM baseline AS static_build
RUN cd ./grackle \
    && cmake -Bbuild -GNinja -DBUILD_SHARED_LIBS=OFF \
    && cmake --build build


# install the static lib
FROM static_build AS static_install
RUN mkdir -p $LOCAL_LIB \
    && cd ./grackle \
    && cmake --install ./build --prefix ${LOCAL_LIB} \
    && rm -r build

# install the static lib variant b
# (specify installation prefix at configure time)
FROM baseline AS static_install_conftimeprefix
RUN mkdir -p $LOCAL_LIB \
    && cd ./grackle \
    && cmake -Bbuild -GNinja -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=$LOCAL_LIB \
    && cmake --build build \
    && cmake --install ./build \
    && rm -r build

# strictly speaking, there's another variant of installation that uses the
# DESTDIR environment variable at install time, but that gets pretty messy
# -> I believe that DESTDIR is simply pretended to partial paths assumed at
#    configure-time (whereas --prefix overwrites the paths)
# -> this is a task for another time....

# build and install the shared lib after the static lib
# -> important: this exists so we can explicitly confirm that the
#    order of operations is unimportant
FROM static_install AS static_shared_install
RUN cd ./grackle \
    && cmake -Bbuild-static -GNinja -DBUILD_SHARED_LIBS=OFF \
    && cmake --build build-static \
    && cmake --install ./build-static --prefix ${LOCAL_LIB}
