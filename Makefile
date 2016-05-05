all: codes build serve

serve: build

.PHONY: codes build serve update

codes:
	@echo codes
	cd codes && make
build:
	jekyll build

watch:
	jekyll build --watch

serve:
	jekyll serve


.PHONY: notebook

notebook:
	rsync -av --exclude='Untitle*' /Users/junkoda/Research/github/mockgallib/notebook/*.html codes/mockgallib/notebook/

