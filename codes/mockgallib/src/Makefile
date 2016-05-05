DIR  := /Users/junkoda/Research/github/mockgallib/libs
SRC  := $(notdir $(wildcard $(DIR)/*.cpp))
SRC  += $(notdir $(wildcard $(DIR)/*.h))

MD   := $(addsuffix .md, $(basename $(SRC)))

.PHONY: all
all:
	cd _src && make
	make $(MD)


%.md: _src/template.erb _src/%.h _src/%.cpp 
	cat $< | sed "s/<%= name %>/$*/" > $@

%.md: _src/template_h.erb _src/%.h 
	cat $< | sed "s/<%= name %>/$*/" > $@

%.md: _src/template_cpp.erb _src/%.cpp 
	cat $< | sed "s/<%= name %>/$*/" > $@
