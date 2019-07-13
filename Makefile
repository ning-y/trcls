github-pages:
	sphinx-build -M html sphinx sphinx-docs
	rm -rf docs && mkdir docs && mv sphinx-docs/html/* docs/
	echo "trcls.ningyuan.io" > docs/CNAME
	touch docs/.nojekyll
