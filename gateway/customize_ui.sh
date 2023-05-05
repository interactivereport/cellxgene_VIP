# Example of how to customize template style: make the index header bright green
# from https://github.com/Novartis/cellxgene-gateway/tree/79_add_docker_example/examples/customized_docker_image
#
# run "CELLXGENE_GATEWAY_DIR=</your/env/lib/python/dir>/site-packages/cellxgene_gateway . ./customize_ui.sh"
#
find "${CELLXGENE_GATEWAY_DIR}/templates" -name index.html -exec sed -i -e 's/<head>/<head>\
> <style> header h3 {color: #0F0;} <\/style>/g' {} \;
