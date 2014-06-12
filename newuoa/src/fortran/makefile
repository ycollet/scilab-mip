main: bigden.o biglag.o calfun.o main.o newuob.o trsapp.o update.o newuoa.o
	g77 -o main bigden.o biglag.o calfun.o main.o newuob.o trsapp.o update.o  newuoa.o

bigden.o: bigden.f
	g77 -g  -c bigden.f
biglag.o: biglag.f
	g77 -g  -c biglag.f
calfun.o: calfun.f
	g77 -g  -c calfun.f
main.o: main.f
	g77 -g  -c main.f
trsapp.o: trsapp.f
	g77 -g  -c trsapp.f
newuob.o: newuob.f
	g77 -g  -c newuob.f
update.o: update.f
	g77 -g  -c update.f
newuoa.o: newuoa.f
	g77 -g  -c newuoa.f