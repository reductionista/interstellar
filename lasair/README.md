
# Lasair client installation

pip3 install lasair

# Lasair client authentication 

For authentication, create ../.venv/lib/python3.12/site-packages/settings.py
 with a single line, setting the string:

lasair_token = [[ LASAIR_API_KEY ]]

api key can be obtained from lasair user profile after creating an account

# Lasair tutorial notebooks

git clone https://github.com/lsst-uk/lasair-examples.git

# Lasair LSST documentation for using the notebooks

https://lasair-lsst.readthedocs.io/en/main/core_functions/python-notebooks.html

# Starting the client

The client is a thin wrapper around confluent_kafka.Consumer. The Kafka server is lasair-lsst-kafka.lsst.ac.uk:9092 and the topic for
receiving all unlinked events is lasair_965nDiaSources5.

To start consuming:

  source .venv/bin/activate
  python python_lib/lasair_kafka_subscriber.py

Or, if direnv is active simply:
  python -m lasair_kafka_subscriber
