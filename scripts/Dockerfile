FROM python:3.10-slim-bullseye

RUN apt update && apt install -y gcc curl musl-dev g++ libhdf5-dev gfortran libffi-dev libssl-dev libopenblas-dev \
	python3-matplotlib python3-scipy python3-numpy python3-pandas python3-numpy libpcre++-dev 

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && . "$HOME/.cargo/env"

WORKDIR /app

COPY . /app/scripts/

RUN mkdir -p /app/packages && chmod u+x /app/scripts/build-deps.sh

VOLUME [ "/app/build" ]

CMD [ "/app/scripts/build-deps.sh" ]