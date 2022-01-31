# Stage 0 - Create from Python3.9.5 image
FROM python:3.9.7-slim-buster as stage0

# Stage 2 - Create virtual environment and install dependencies
FROM stage0 as stage1
COPY requirements.txt /app/requirements.txt
RUN /usr/local/bin/python3 -m venv /app/env
RUN /app/env/bin/pip install -r /app/requirements.txt

# Stage 1 - Copy MetroMan code
FROM stage1 as stage2
COPY ./moi /app/moi/

# Stage 3 - Execute algorithm
FROM stage2 as stage3
COPY run_MOI.py /app/run_MOI.py
LABEL version="1.0" \
	description="Containerized MOI algorithm." \
	"confluence.contact"="ntebaldi@umass.edu" \
	"algorithm.contact"="durand.8@osu.edu"
ENTRYPOINT ["/app/env/bin/python3", "/app/run_MOI.py"]