FROM brinkmanlab/psortb_commandline:1.0.2

RUN apt-get update && apt-get install -y \
    'python3' \
    'python3-pip' \
    'hmmer' \
    'prodigal'

RUN pip3 install --upgrade \
    'pip' \
    'biopython' \
    'networkx==2.2' \
    'conda' \
    'Keras==0.3.3' \
    'tensorflow==1.1.0' \
    'ptvsd' \
    'h5py'

# RUN wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && \
    # gunzip Pfam-A.hmm.gz

ADD . ./

# Need to have this downloaded beforehand if building image locally - they send an email to the source.
# RUN tar -xvzf hmmto   p_2.1.tgz && rm hmmtop_2.1.tgz
# RUN cd hmmtop_2.1 && cc hmmtop.c -lm -o hmmtop
# RUN hmmpress db/Pfam-A.hmm
ADD . ./

ENTRYPOINT []