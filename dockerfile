FROM debian:bullseye
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update 
RUN apt-get install -y build-essential protobuf-compiler libzmq3-dev curl g++  -y git
RUN apt install -y libprotobuf-dev protobuf-compiler
RUN apt-get install -y software-properties-common gcc #&&  add-apt-repository -y ppa:deadsnakes/ppa
RUN apt-get update && apt-get install -y python3-pip swig libfftw3-dev libfftw3-doc libssl-dev cmake 
RUN export CUDA_PATH=/usr
RUN pip install numpy matplotlib pandas plotly 
RUN git clone https://github.com/WinklerTB/MicMag2
RUN mkdir MicMag2/src/build
WORKDIR MicMag2/src/build
RUN cmake .. 
RUN make -j 4

