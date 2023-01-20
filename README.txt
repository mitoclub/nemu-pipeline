# Download refseq from 
https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/

# extract gorilla from refseq file
```
awk 'BEGIN {RS=">"} /Gorilla/ {printf ">"$o}' mitochondrion.1.1.genomic.fna > gorilla.fna
```

# if-else in bash
if [ `grep -c korilla hello.pl` -gt 1 ]; then echo yes; else echo no; fi

# create nucleotide blast-db
makeblastdb -in humans_mt.fasta -dbtype nucl -title "Human mtDNA"
sudo singularity exec --bind /mnt/data/export:/export src/image_pipeline-2.5.sif makeblastdb -in /export/data/humans_mt.fasta -dbtype nucl -title "Human mtDNA" -parse_seqids

# download genes to create db
mkdir DIR
perl /opt/scripts/geneUploader.pl SPECIES GENE DIR

# change ownership of dir
sudo chown -R --reference=data/gorilla.fna data/base

# mount data-disk
sudo mount /dev/mapper/vg-data /mnt/data



when:
hits_file.toString() != "hits_yes.txt"

validExitStatus 1
errorStrategy 'finish'

# KG query to NCBI
mtDNA vertebrate


# run container
docker run --privileged -m 10G -p 8080:80 -v /mnt/data/export:/export -dti ummsbiocore/dolphinnext-studio
#docker run --storage-opts dm.basesize=50G -m 10G -p 8080:80 -v ~/export:/export -dti ummsbiocore/dolphinnext-studio

# startup dolphinnext
docker exec -it $(docker ps | tail -n 1 | cut -d " " -f 1) bash

# change default hostname of web app in config file  BASE_PATH & PUBWEB_URL
vim /export/dolphinnext/config/.sec

# change ssl configs
vim /etc/apache2/sites-enabled/default-ssl.conf
                SSLCertificateFile      /etc/ssl/certs/ssl-cert-snakeoil.pem
                SSLCertificateKeyFile /etc/ssl/private/ssl-cert-snakeoil.key

                SSLCertificateFile      /export/pki/kantiana_ru_2023_12_10.crt
                SSLCertificateKeyFile /export/pki/kantiana_ru_2023_12_10.key

sudo service apache2 restart

# our host
bioinfo.int.kantiana.ru

# change user for /export dir inside docker to run pipelines
chown -v docker /export


```
`tail -n $(($(wc -l cox2.faa | cut -f 1 -d " ") - 1))`

# make db
## DB must contain more than 2 species!!!

mkdir DIR
perl /opt/scripts/geneUploader.pl SPECIES GENE DIR
```


# 16.4, page 71
• Optional Inputs:
If you want use optional input parameter, you can check the optional checkbox. This feature allows flexibility for the
user while defining process since the process will be executed in spite of the absence of the input parameter. Please
check the example below for the use case:

# directives to processes
https://www.nextflow.io/docs/latest/process.html?highlight=clause#directives
useful:
   errorStrategy
   time
   validExitStatus
   queue
   pod


# Когда запускал локально на mr@Smith

docker run --privileged -m 5G -p 8080:80 -v ~/dolphin/exportm:/export -dti ummsbiocore/dolphinnext-studio
docker exec -it $(docker ps | tail -n 1 | cut -d " " -f 1) bash

# далее внутри контейнера
startup  # не создался нормальный ключ для локального запуска в контейнере, точнее создался, вроде, но в приложении были false
   15  cd /home/docker/
   16  cd .ssh/
   17  ssh-keygen -t rsa
   18  ls
   19  cat id_rsa  # copied to private in key creation menu
   20  cat id_rsa.pub  # copied to public (после этого назначил в дефолтном профиле этот ключ, fastq-файлы стали видны)
   21  #chown -v docker /export/
   22  ll /export/
   23  chown -v docker /export/  # не хотело запускаться, пришлось изменить пользователя, группа не изменилась (root)
   24  ll /export/
   25  history

# containers
1. p.def - only python and its programms
2. pipeline-2.5.def - working through p.def result container
3. pipeline.def - independent and working

notes:
- Processes should have different names
- one process run one time
- types of processes outputs must be different
- 