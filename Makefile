all: build serve

serve: build

.PHONEY: build serve update

build:
	bundle exec jekyll build

watch:
	bundle exec jekyll build --watch

serve:
	bundle exec jekyll serve

update:
	bundle update


.PHONY: notebook

notebook:
	rsync -av --exclude='Untitle*' /Users/junkoda/Research/github/mockgallib/notebook/*.html codes/mockgallib/notebook/

