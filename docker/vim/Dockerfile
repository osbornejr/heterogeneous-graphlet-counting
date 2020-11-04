FROM julia:latest

USER root
RUN apt-get update && apt-get install -y \
	vim \
	git 
COPY .vim $HOME/.vim/ 
USER $USER
ENV TERM xterm-256color
#CMD ["vim", "test"]
CMD ["/bin/bash"]


