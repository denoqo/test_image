
build:
	DOCKER_BUILDKIT=1 docker build --no-cache -t eas-test .

run:
	docker run -it --rm -p 3838:3838 eas-test

runmap:
	sudo docker container run --rm  -p 3838:3838 \
	-v //home/ubuntu/test_image/shiny-app/app/:/srv/shiny-server/test_app \
	-v /home/ubuntu/log/shiny-server/:/var/log/shiny-server/ \
	qx658/eas-shiny
