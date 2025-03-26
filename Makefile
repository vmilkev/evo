.PHONY : all
all:
	cd libplinkio && make
	cd evoped && make
	cd evolm && make
	cd evogen && make

.PHONY : test
test:
	cd evoped && make $@
	cd evolm && make $@
	cd evogen && make $@

.PHONY : rclean
rclean:
	cd evolm && make $@
	cd evoped && make $@
	cd evogen && make $@

.PHONY : clean
clean:
	cd libplinkio && make $@
	cd evoped && make $@
	cd evolm && make $@
	cd evogen && make $@