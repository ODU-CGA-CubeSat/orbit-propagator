#!/usr/bin/env bash

asciidoctor research/OrbitPropagation.adoc -o dist/OrbitPropagation.pdf -r asciidoctor-pdf -r asciidoctor-diagram -r asciidoctor-bibtex -b pdf

asciidoctor research/OrbitPropagation.adoc -o dist/OrbitPropagation.html -r asciidoctor-diagram -r asciidoctor-bibtex