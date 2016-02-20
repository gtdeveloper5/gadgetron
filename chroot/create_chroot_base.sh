#!/bin/bash

if [ $(id -u) -ne 0 ]; then
    echo -e "\nPlease start the script as a root or sudo!\n"
    exit 1
else
    if [ $# -ge 1 ]; then

        # --ARGUMENTS-- (example)

        # CHROOT_GADGETRON_BINARY_DIR:        /home/ubuntu/gadgetron/build

        # -----------------------------------------------------------------------------------
        # input parameters
        # -----------------------------------------------------------------------------------

        CHROOT_GADGETRON_BINARY_DIR=${1}
        echo CHROOT_GADGETRON_BINARY_DIR: ${CHROOT_GADGETRON_BINARY_DIR}

        rm -rf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
        mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
        mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups

        # install debian bootstrap
        apt-get install debootstrap -y
        debootstrap --variant=buildd --arch amd64 trusty ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron http://gb.archive.ubuntu.com/ubuntu/

        # install python libraries
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get install --no-install-recommends --no-install-suggests --yes wget build-essential emacs python-pip libhdf5-serial-dev cmake git-core libboost-all-dev libfftw3-dev h5utils hdf5-tools liblapack-dev libxml2-dev libfreetype6-dev pkg-config libxslt-dev libarmadillo-dev libace-dev gcc-multilib libgtest-dev liblapack-dev libatlas-base-dev libatlas-dev libplplot-dev supervisor
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get install software-properties-common python-libxml2 -y
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron pip install Cython
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron pip install numpy pyxb h5py psutil Twisted lxml matplotlib  

        TAR_FILE_NAME=gadgetron-base-`date '+%Y%m%d-%H%M'`

        tar -zcf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz" --directory "${CHROOT_GADGETRON_BINARY_DIR}/chroot" --exclude=./chroot-root/gadgetron/var --exclude=./chroot-root/gadgetron/dev --exclude=./chroot-root/gadgetron/sys --exclude=./chroot-root/gadgetron/proc --exclude=./chroot-root/gadgetron/root ./chroot-root

        rm -rf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root"

        chmod 666 "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz"
        exit 0
    else
        echo -e "\nUsage:  $0 (gadgetron binary dir)\n"
        exit 1
    fi
fi
