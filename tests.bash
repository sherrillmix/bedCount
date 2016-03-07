#!/bin/bash

set -e

./bedCount 2>/dev/null && { echo "Missing files did not fail"; exit 1; }
./bedCount -h |grep Usage|| { echo "-h did not generate usage"; exit 1; }
./bedCount -Z && { echo "Weird arg did not fail"; exit 1; }
./bedCount -$'\05' && { echo "Weird arg did not fail"; exit 1; }

./bam2depth 2>/dev/null && { echo "Missing files did not fail"; exit 1; }
./bam2depth -h |grep Usage|| { echo "-h did not generate usage"; exit 1; }
./bam2depth -Z && { echo "Weird arg did not fail"; exit 1; }
./bam2depth -$'\05' && { echo "Weird arg did not fail"; exit 1; }
