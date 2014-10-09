CCTYPE = -DGCC # sinon -DSUN
GCC = gcc -ggdb -Wall #-pg
GPP = g++ -ggdb -Wall #-pg

OBJS = LL_Lolimot.o \
	LL_Mesure.o \
	LL_Cut.o \
	LL_Partition.o \
	LL_Dimension.o \
	DefaultParam.o \
	ExportCommonC_Cpp.o \
	ExportC.o \
	ExportCpp.o \
	ExportMatlab.o \
	ExportMinMaxMap.o \
	Filter.o \
	Learning.o \
	DisplayData.o \
	Analysis.o \
	Post.o \
	Validate.o \
	libExport.a \
	libTrain.a \
	libLolimot.a \
	lolimot

all: $(OBJS) 

LL_Lolimot.o: LL_Lolimot.h LL_Lolimot.cpp
	$(GPP) -c LL_Lolimot.cpp -o LL_Lolimot.o -I. -I./SmallLapack -I./Lasso
LL_Mesure.o: LL_Mesure.h LL_Mesure.cpp
	$(GPP) -c LL_Mesure.cpp -o LL_Mesure.o -I.
LL_Cut.o: LL_Cut.h LL_Cut.cpp
	$(GPP) -c LL_Cut.cpp -o LL_Cut.o -I.
LL_Partition.o: LL_Partition.h LL_Partition.cpp
	$(GPP) -c LL_Partition.cpp -o LL_Partition.o -I.
LL_Dimension.o: LL_Dimension.h LL_Dimension.cpp
	$(GPP) -c LL_Dimension.cpp -o LL_Dimension.o -I.
DefaultParam.o: DefaultParam.cpp
	$(GPP) -c DefaultParam.cpp -o DefaultParam.o -I.
ExportCommonC_Cpp.o: ExportCommonC_Cpp.cpp ExportCommonC_Cpp.h
	$(GPP) -c ExportCommonC_Cpp.cpp -o ExportCommonC_Cpp.o -I.
ExportC.o: ExportC.cpp ExportC.h
	$(GPP) -c ExportC.cpp -o ExportC.o -I.
ExportCpp.o: ExportCpp.cpp ExportCpp.h
	$(GPP) -c ExportCpp.cpp -o ExportCpp.o -I.
ExportMatlab.o: ExportMatlab.cpp ExportMatlab.h
	$(GPP) -c ExportMatlab.cpp -o ExportMatlab.o -I.
ExportMinMaxMap.o: ExportMinMaxMap.cpp ExportMinMaxMap.h
	$(GPP) -c ExportMinMaxMap.cpp -o ExportMinMaxMap.o -I.
Filter.o: Filter.cpp Filter.h
	$(GPP) -c Filter.cpp -o Filter.o -I.
Learning.o: Learning.cpp Learning.h
	$(GPP) -c Learning.cpp -o Learning.o -I.
DisplayData.o: DisplayData.cpp DisplayData.h
	$(GPP) -c DisplayData.cpp -o DisplayData.o -I.
Analysis.o: Analysis.cpp Analysis.h
	$(GPP) -c Analysis.cpp -o Analysis.o -I.
Post.o: Post.cpp Post.h
	$(GPP) -c Post.cpp -o Post.o -I.
Validate.o: Validate.cpp Validate.h
	$(GPP) -c Validate.cpp -o Validate.o -I.
TrainLolimotStruct.o: TrainLolimotStruct.cpp TrainLolimotStruct.h
	$(GPP) -c TrainLolimotStruct.cpp -o TrainLolimotStruct.o -I.
libExport.a: ExportCommonC_Cpp.o ExportC.o ExportCpp.o ExportMatlab.o
	ar -cru libExport.a ExportCommonC_Cpp.o ExportC.o ExportCpp.o ExportMatlab.o
libTrain.a: DisplayData.o Learning.o Filter.o Analysis.o Validate.o TrainLolimotStruct.o
	ar -cru libTrain.a DisplayData.o Learning.o Filter.o Analysis.o Validate.o TrainLolimotStruct.o
libLolimot.a: LL_Lolimot.o LL_Partition.o LL_Mesure.o LL_Dimension.o DefaultParam.o LL_Cut.o Post.o
	ar -cru libLolimot.a LL_Lolimot.o LL_Partition.o LL_Mesure.o LL_Dimension.o DefaultParam.o LL_Cut.o Post.o
lolimot: TrainLolimot2.cpp libLolimot.a libTrain.a libExport.a
	cd Lasso; make all
	cd SmallLapack; make all
	$(GPP) TrainLolimot2.cpp -o lolimot -I. libLolimot.a ./Lasso/lasso.o libTrain.a libExport.a -lm ./SmallLapack/liblapack.a $(CCTYPE)
clean:
	rm *.o
	cd SmallLapack; make clean
	cd Lasso; make clean
