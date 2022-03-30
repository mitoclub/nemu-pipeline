# create nucleotide blast-db
makeblastdb -in humans_mt.fasta -dbtype nucl -title "Human mtDNA"
sudo singularity exec --bind /mnt/data/export:/export src/image_pipeline-2.5.sif makeblastdb -in /export/data/humans_mt.fasta -dbtype nucl -title "Human mtDNA"


# mount data-disk
sudo mount /dev/mapper/vg-data /mnt/data


# run container
docker run --privileged -m 10G -p 8080:80 -v /mnt/data/export:/export -dti ummsbiocore/dolphinnext-studio
#docker run --storage-opts dm.basesize=50G -m 10G -p 8080:80 -v ~/export:/export -dti ummsbiocore/dolphinnext-studio

# startup dolphinnext
docker exec -it $(docker ps | tail -n 1 | cut -d " " -f 1) bash

# change default hostname of web app in config file  BASE_PATH & PUBWEB_URL
vim /export/dolphinnext/config/.sec

# our host
bioinfo.int.kantiana.ru

# change user for /export dir inside docker to run pipelines
chown -v docker /export



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

