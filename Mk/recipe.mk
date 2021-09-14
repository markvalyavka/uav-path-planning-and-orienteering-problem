uniq = $(if $1,$(firstword $1) $(call uniq,$(filter-out $(firstword $1),$1)))

CXXFLAGS:=$(call uniq,$(CXXFLAGS))

clean_make: clean $(TARGET)

bin: $(TARGET)

$(OBJS): %.o: src/%.cpp
	echo "making $<"
	$(CXX) -c $< $(CXXFLAGS) $(CPPFLAGS) -o $(OBJ_DIR)/$@

$(TARGET): create_directories  obj_subdirs $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ_DIR)/*.o $(LDFLAGS)


all_subdirs:    
	echo $(OPSYS)
	@for i in $(SUBDIRS) ;\
	do \
	echo "making" all "in $(CURRENT_DIR)/$$i..."; \
	$(MAKE) -C $$i all || exit 1 ;\
	done

obj_subdirs: all_subdirs
	echo "Copy objs"
	@for i in $(SUBDIRS) ;\
	do \
	$(MAKE) --ignore-errors -C $$i obj ;\
	echo "coping all in $(CURRENT_DIR)/$$i..."; \
	cp $$i/*.o $(OBJ_DIR)/; \
	done

create_directories:
	echo "create dircetory $(OBJ_DIR)"
	mkdir -p $(OBJ_DIR)

dependencies: agilicious_lib

agilicious_lib:
	echo "making agilicious"
	cd ./agilicious/agilib/build; cmake -DCMAKE_CXX_FLAGS="-Werror" -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TEST=ON -DEIGEN_FROM_SYSTEM=OFF .. ; make -j8;
	#cd ./agilicious/agilib/build; cmake .. ; make -j8; #-DUNSAFE_MATH=OFF -DCMAKE_VERBOSE_MAKEFILE=ON
	#cd ./agilicious/agilib/build; cmake .. -DCMAKE_BUILD_TYPE=Debug; make;
	#-DCMAKE_BUILD_TYPE=Release

clean_dependencies: clean_agilicious_lib

clean_agilicious_lib:
	echo "cleaning flann"
	rm -rf ./lib/agilicious/agilib/build/*

clean:
	$(RM) $(TARGET)
	$(RM) $(OBJ_DIR)/*
	@for i in $(SUBDIRS) ;\
        do \
        echo "cleaning" all "in $(CURRENT_DIR)/$$i..."; \
        $(MAKE) -C $$i clean; \
        done

clean_all: clean clean_dependencies

