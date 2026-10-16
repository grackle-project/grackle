# we refer to the various "layers" in this file as "image-targets"
# - each image-target exists to be used to create a container for running an
#   installtest test-case.
# - we encode a basic summary in comments just before the declaration of each
#   "image-target". Such a description is enclosed between a line that says
#   `# DESCRIPTION-START` and a line that says `# DESCRIPTION-END`

# currently, it's important to invoke docker from the root directory in order to copy
# over the source directory

# DESCRIPTION-START
# A basic image where the Grackle source directory has been copied and dependencies have
# been installed. No builds have been configured yet
# DESCRIPTION-END
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

COPY --chown=$username ./ ./grackle

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

# DESCRIPTION-START
# The shared library form of Grackle is built and can be found in its build directory,
# ~/grackle/build
# DESCRIPTION-END
FROM baseline AS shared_build
RUN cd ./grackle \
    && cmake -Bbuild -GNinja -DBUILD_SHARED_LIBS=ON \
    && cmake --build build


# DESCRIPTION-START
# The shared library form of Grackle has been installed (and the build directory was
# cleaned up)
# DESCRIPTION-END
FROM shared_build AS shared_install
RUN mkdir -p $LOCAL_LIB \
    && cd ./grackle \
    && cmake --install ./build --prefix ${LOCAL_LIB} \
    && rm -r build


# DESCRIPTION-START
# Both the shared and shared library forms of Grackle have been built and installed
# (the shared library was installed first).
#
# Part of the reason this exists is so that we can explicitly confirm that the
# order of operations is unimportant
# DESCRIPTION-END
FROM shared_install AS shared_static_install
RUN cd ./grackle \
    && cmake -Bbuild-static -GNinja -DBUILD_SHARED_LIBS=OFF \
    && cmake --build build-static \
    && cmake --install ./build-static --prefix ${LOCAL_LIB} \
    && rm -r ./build-static


# DESCRIPTION-START
# The static library form of Grackle is built and can be found in its build directory,
# ~/grackle/build
# DESCRIPTION-END
FROM baseline AS static_build
RUN cd ./grackle \
    && cmake -Bbuild -GNinja -DBUILD_SHARED_LIBS=OFF \
    && cmake --build build


# DESCRIPTION-START
# The static library form of Grackle has been installed (and the build directory was
# cleaned up).
# DESCRIPTION-END
FROM static_build AS static_install
RUN mkdir -p $LOCAL_LIB \
    && cd ./grackle \
    && cmake --install ./build --prefix ${LOCAL_LIB} \
    && rm -r build

# DESCRIPTION-START
# The static library form of Grackle has been installed (and the build directory was
# cleaned up). In contrast to all other image-targets that install Grackle, this
# image-target sets the installation directory when the build is initially configured by
# setting the standard ``CMAKE_INSTALL_PREFIX`` cmake variable.
#
# For added context, all of the other image-targets specify the desired installation
# direction when executing the ``cmake`` program by passing the path through via the
# ``--prefix`` flag.
# DESCRIPTION-END
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

# DESCRIPTION-START
# Both the shared and shared library forms of Grackle have been built and installed
# (the static library was installed first).
#
# Part of the reason this exists is so that we can explicitly confirm that the
# order of operations is unimportant
# DESCRIPTION-END
FROM static_install AS static_shared_install
RUN cd ./grackle \
    && cmake -Bbuild-shared -GNinja -DBUILD_SHARED_LIBS=ON \
    && cmake --build build-shared \
    && cmake --install ./build-shared --prefix ${LOCAL_LIB}
