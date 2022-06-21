all: ProfileMakers Archive

ProfileMakers::
	(cd ProfileMakers; make)
	(mvexist.sh ProfileMakers/src/ProfileMakers_Dict_rdict.pcm lib/)
	(mvexist.sh ProfileMakers/libProfileMakers.rootmap lib/)

Archive::
	(tar -zcf lib/DataFormats.tar.gz ProfileMakers)

clean::
	(cd ProfileMakers; make clean)

distclean::
	(cd ProfileMakers; make distclean)

LibTarFile = tar/lib.tar.gz
$(LibTarFile): $(wildcard ./lib/*)
	tar -czf $@ ./lib/*

DataTarFile = tar/data.tar.gz
DataFiles = $(shell find data/$(LArProfV)/ -type f -name '*')
$(DataTarFile): $(DataFiles)
	tar --exclude=data/$(LArProfV)/Sample -czf $@ ./data/$(LArProfV)/

CondorTar: $(LibTarFile) $(DataTarFile)
