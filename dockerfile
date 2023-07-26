FROM debian:bullseye
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update 
RUN apt-get install -y build-essential protobuf-compiler libzmq3-dev curl g++  -y git
RUN apt install -y libprotobuf-dev protobuf-compiler
RUN apt-get install -y software-properties-common gcc #&&  add-apt-repository -y ppa:deadsnakes/ppa
RUN apt-get update && apt-get install -y python3-pip swig fftw3 libssl-dev cmake 
RUN pip install numpy matplotlib pandas plotly 
RUN mkdir MM2
RUN git clone https://github.com/WinklerTB/MicMag2 
#WORKDIR src/build_orig
#RUN cmake ..
#RUN make -j 4

