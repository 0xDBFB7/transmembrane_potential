wget --continue www.coin-or.org/download/source/ADOL-C/ADOL-C-2.6.3.tgz
tar zxvf ADOL-C-2.6.3.tgz
cd ADOL-C-2.6.3
mkdir ./ThirdParty
cd ./ThirdParty
wget --continue http://archive.ubuntu.com/ubuntu/pool/universe/c/colpack/colpack_1.0.10.orig.tar.gz
tar zxvf colpack_1.0.10.orig.tar.gz
mv ColPack-1.0.10 ColPack
cd ColPack
./configure --prefix=/usr/local
make
sudo make install
cd $HOME/ADOL-C-2.6.3
./configure --prefix=/usr/local --enable-sparse --with-colpack=/usr/local
make
sudo make install
