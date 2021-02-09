FROM gitlab.c3s.unito.it:5000/gbirolo/ngs_basic:v1.2 as build
RUN apt-get update && apt-get install -y --no-install-recommends unzip wget

#COPY Binaries/gatk-4.1.9.0.zip /
#RUN unzip gatk-4.1.9.0.zip
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip -O gatk.zip && unzip gatk.zip

#COPY Binaries/snpEff_latest_core.zip /
#RUN unzip snpEff_latest_core.zip
RUN wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip/download -O snpeff.zip && unzip snpeff.zip

RUN wget https://github.com/brentp/vcfanno/releases/download/v0.3.2/vcfanno_linux64

FROM gitlab.c3s.unito.it:5000/gbirolo/ngs_basic:v1.2

COPY --from=build /gatk-4.1.9.0 /opt/gatk
RUN ln -s /opt/gatk/gatk /opt/bin/

COPY --from=build /snpEff /opt/snpEff
# change snpEff default data directory for occam
RUN sed -i 's@^data.dir = ./data/$@data.dir = /scratch/shared/ngs_geno/snpEff_data@' /opt/snpEff/snpEff.config

COPY --from=build /vcfanno_linux64 /opt/bin/vcfanno
#COPY Binaries/vcfanno_linux64 /opt/bin/vcfanno
RUN chmod a+x /opt/bin/vcfanno

COPY Pipelines/* /opt/pipelines/
