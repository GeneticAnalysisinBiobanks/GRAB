# Build the Docker image
sudo docker build -t geneticanalysisinbiobanks/grab:v0.2.3 - <<EOF
FROM condaforge/miniforge3
RUN conda install r-grab
EOF

# Upload the image to Docker Hub
sudo docker login -u geneticanalysisinbiobanks
sudo docker push geneticanalysisinbiobanks/grab:v0.2.3

# Test the image
sudo docker pull geneticanalysisinbiobanks/grab:v0.2.3
sudo docker run --name=miaolin-test1 --cpus=1 --memory=1G --rm geneticanalysisinbiobanks/grab:v0.2.3 R -e "library(GRAB); message('GRAB loaded successfully')"

# Clean up
sudo docker rmi 307a9583dfea
sudo docker system prune -a --volumes
