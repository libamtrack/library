#!/bin/bash
if [[ "$TRAVIS_PULL_REQUEST" == "false" && "$TRAVIS_BRANCH" == "master" ]]
then
  curl -s -X POST -H "Content-Type: application/json" -H "Accept: application/json" -H "Travis-API-Version: 3" -H "Authorization: token $TRAVIS_TOKEN" -d '{ "request": { "branch":"master", "message": "Build Travis after changes in libamtrack/library" }}' https://api.travis-ci.com/repo/libamtrack%2FDockerFiles/requests
fi
