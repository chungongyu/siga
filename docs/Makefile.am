if HELP2MAN
man_MANS=siga.1
siga.1:
	help2man -o siga.1 ${top_builddir}/src/siga${EXEEXT}
else
man_MANS=
endif !HELP2MAN

SUBDIRS=_layouts
EXTRA_DIST=404.html \
        Gemfile \
        _config.yml \
        favicon.png \
        index.md \
        logo_with_text.png
html-env:
	bundle install --path=vendor/bundle
html-local: html-env
	bundle exec jekyll build

clean-local:
	rm -rf _site siga.1
distclean-local:
	rm -rf vendor
