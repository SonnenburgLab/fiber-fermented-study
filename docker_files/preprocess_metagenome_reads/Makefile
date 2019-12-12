NAME=preprocess_metagenome_reads
TAG:=$(shell /bin/date +%Y%m%d%H%M%S)
REGISTRY=hwastyk/$(NAME)

all: build push run

build:
	docker image build -t $(NAME):$(TAG) -t $(NAME):latest -f Dockerfile . # add --no-cache after build if needed
	docker image tag $(NAME):latest $(REGISTRY):latest
	docker image tag $(NAME):$(TAG) $(REGISTRY):$(TAG)
	
push:	
	docker image push $(REGISTRY):$(TAG)
	docker image push $(REGISTRY):latest
	echo "bmerrill9/$(NAME):$(TAG)" > LATEST

run:
	docker container run --rm $(NAME):latest
