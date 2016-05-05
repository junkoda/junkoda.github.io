
DIR  := /Users/junkoda/Research/github/mockgallib/libs
SRC  := $(notdir $(wildcard $(DIR)/*.cpp))
SRC  += $(notdir $(wildcard $(DIR)/*.h))
HTML := $(addsuffix .html, $(basename $(SRC)))

all: $(SRC)

%.cpp: $(DIR)/%.cpp
	cp $< $@

%.h: $(DIR)/%.h
	cp $< $@




