
#.PHONY: all sim chisel clean
.PHONY: all sim doc clean

#all: chisel sim
all: sim

#sim: chisel
sim:
	cd /home/kiviukko/digitaltransmitter/thesdk_template/Entities/Transmitter/Transmitter && \
	python3 __init__.py

#chisel:
#	cd /home/kiviukko/digitaltransmitter/thesdk_template/Entities/Transmitter/chisel && \
#	make

doc:
	cd /home/kiviukko/digitaltransmitter/thesdk_template/Entities/Transmitter/doc && \
	make html

clean:
	#cd /home/kiviukko/digitaltransmitter/thesdk_template/Entities/Transmitter/chisel && \
	#make clean && \
	rm -rf /home/kiviukko/digitaltransmitter/thesdk_template/Entities/Transmitter/simulations/* 
   
