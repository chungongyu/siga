man_MANS=siga.1
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
	rm -rf _site
distclean-local:
	rm -rf vendor
