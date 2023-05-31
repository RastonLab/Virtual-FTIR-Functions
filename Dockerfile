FROM python:3.10-slim-bullseye as base

RUN apt update && apt install wget tar -y

WORKDIR /app

COPY scripts /app/scripts/

RUN mkdir -p /app/packages && wget https://objectstorage.us-phoenix-1.oraclecloud.com/n/ax3qqxs2wjil/b/prebuilt-deps/o/virtual-ftir-functions-deps-$(uname -m)-py3.10.tar.gz

RUN tar -xvf virtual-ftir-functions-deps-$(uname -m)-py3.10.tar.gz -C /app/packages

RUN pip3 install --no-index --find-links=/app/packages -r /app/scripts/requirements.txt

FROM python:3.10-slim-bullseye

# copy pip packages from base
COPY --from=base /usr/local/lib/python3.10/site-packages /usr/local/lib/python3.10/site-packages

COPY ./*.py .

CMD ["python", "app.py"]
