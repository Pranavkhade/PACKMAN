FROM centos:7

# System dependencies
RUN yum -y install git which python-virtualenv gcc gcc-c++

# SSH deployment key
RUN mkdir ~/.ssh
COPY ppackman-key ~/.ssh/
RUN chmod 400 ~/.ssh/ppackman-key

# Setup virtualenv
RUN virtualenv ppackman-env
RUN source ppackman-env/bin/activate

# Clone the source
RUN git clone git@github.com:ResearchIT/PPACKMAN.git

# Install python requirements
WORKDIR PPACKMAN
RUN pip install -r requirements.txt

