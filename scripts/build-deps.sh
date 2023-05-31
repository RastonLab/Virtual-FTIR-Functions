#!/bin/bash
# This script is used to build the dependencies for the project.

set -e
# print usage
usage() {
	echo "Usage: $0 [--arm64]" 1>&2
	echo "Builds the dependencies for the project." 1>&2
	echo "  --help, -h:" 1>&2
	echo "    Prints this help message." 1>&2
	echo "  --arm64:" 1>&2
	echo "    Builds the dependencies for the project. If --arm64 is passed as an argument, the dependencies will be built for arm64." 1>&2
	echo "    ONLY USE THIS IF YOU ARE NOT BUILDING ON AN ARM MACHINE" 1>&2
	echo "" 1>&2
	echo "ENVIRONMENT VARIABLES:" 1>&2
	echo "  S3_BUCKET:" 1>&2
	echo "    The name of the S3 bucket to upload the dependencies to." 1>&2
	echo "    AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY must be set in order to upload to S3." 1>&2
	echo "  OCI_OBJECT_STORAGE_URL:" 1>&2
	echo "    The URL of the OCI bucket to upload the dependencies to." 1>&2
	echo "    OCI_OBJECT_STORAGE_URL must be pre-authenticated to write to the bucket." 1>&2
	exit 1
}
# print usage if --help or -h is passed as an argument
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
	usage
fi

if [ -f /.dockerenv ]; then
    echo "Containerized environment detected";
else
    echo "Host environment detected"
	# detect if docker is installed
	if ! [ -x "$(command -v docker)" ]; then
		echo "Docker is not installed. Please install docker to continue."
		exit 1
	fi
	IMAGE_NAME="ghcr.io/rastonlab/build-deps"

	# detect if platform build is requested
	if [ "$1" == "--arm64" ]; then
		echo "Building for arm64"
		docker buildx build --platform linux/arm64 -t $IMAGE_NAME --load .
	else
		echo "Single-arch build requested"
		docker build -t $IMAGE_NAME .
	fi
	# if uploading to OCI or S3, don't mount the build directory
	if [ -z "$S3_BUCKET" ] && [ -z "$OCI_OBJECT_STORAGE_URL" ]; then
		BUILD_VOLUME="-v $(pwd)/out:/app/build"
	fi

	docker run --rm $BUILD_VOLUME -e S3_BUCKET=$S3_BUCKET -e OCI_OBJECT_STORAGE_URL=$OCI_OBJECT_STORAGE_URL $IMAGE_NAME
	# clean up
	docker image rm $IMAGE_NAME
	exit 0
fi
# add the cargo bin to the path
source "$HOME/.cargo/env"

# install the python dependencies
pip3 wheel --wheel-dir /app/packages --prefer-binary -r /app/scripts/requirements.txt

# create a compressed tarball of the directory with the name of the architecture and python version
TARBALL_NAME="virtual-ftir-functions-deps-$(uname -m)-py$(python3 -c 'import sys; print(sys.version_info.major,sys.version_info.minor,sep=".")').tar.gz"
cd /app/packages
tar -czvf /app/build/$TARBALL_NAME ./*

# if S3_BUCKET environment variable is set, upload the dependencies to the specified s3 bucket
if [ -z "$S3_BUCKET" ]; then
	echo "S3_BUCKET environment variable not set. Skipping upload."
else
	echo "Uploading dependencies to S3 bucket $S3_BUCKET"
	# install the aws cli
	pip3 install awscli
	# upload the tarball to the s3 bucket
	aws configure set aws_access_key_id $AWS_ACCESS_KEY_ID
	aws configure set aws_secret_access_key $AWS_SECRET_ACCESS_KEY
	aws s3 cp /app/build/$TARBALL_NAME s3://$S3_BUCKET/$TARBALL_NAME 
fi

# if OCI_OBJECT_STORAGE_URL environment variable is set, upload the dependencies to the specified OCI bucket
# NOTE: OCI_OBJECT_STORAGE_URL must be pre-authenticated to write to the bucket
if [ -z "$OCI_OBJECT_STORAGE_URL" ]; then
	echo "OCI_OBJECT_STORAGE_URL environment variable not set. Skipping upload."
else
	echo "Uploading dependencies to OCI bucket"
	curl -X PUT -T /app/build/$TARBALL_NAME $OCI_OBJECT_STORAGE_URL
fi
