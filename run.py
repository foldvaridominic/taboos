import sys
import argparse
import logging
import time
import re

from constants import SPECIAL_CHARACTER_N
from tabooset import TabooSet


def run(organism_regex, enzyme_regex, rsite_regex, max_wildcard):

    idx=0
    start_time = time.time()
    re_organism = re.compile(organism_regex)
    re_enzyme = re.compile(enzyme_regex)
    re_rsite = re.compile(rsite_regex)
    organism = None
    prev_org = None
    data = []
    collect = False
    not_connected_count = 0


    for line in sys.stdin:
        name_match = re_organism.search(line)
        enzyme_match = re_enzyme.search(line)
        rsite_match = re_rsite.search(line)

        if enzyme_match:
            collect = True
        if collect and name_match:
            prev_org = organism if organism is not None else name_match.group('name')
            organism = name_match.group('name')
        if collect and data and prev_org != organism:
            idx += 1
            logger.info("%s name: %s", idx, prev_org)
            ts = TabooSet(data)
            connected = ts.connected
            logger.info("taboo count: %s | %s", len(ts.taboo_strings),
                'not connected' if not connected else 'connected')
            data = []
            prev_org = organism
            if not connected:
                not_connected_count += 1
        if collect and rsite_match:
            restriction_site = rsite_match.group('nucleotides')
            if sum(1 for i in restriction_site if i == SPECIAL_CHARACTER_N) <= max_wildcard:
                current_data = restriction_site
            else:
                current_data = 'too many wildcard characters'
            data.append(current_data)
            collect = False
            prev_org = organism

    if data:
        idx += 1
        logger.info("%s name: %s", idx, prev_org)
        logger.info("last one")
        ts = TabooSet(data)
        connected = ts.connected
        logger.info("taboo count: %s | %s", len(ts.taboo_strings),
            'not connected' if not connected else 'connected')

    logger.info("progress: %d", idx)
    logger.info("finished in %ss", time.time() - start_time)
    logger.info("not connected count: %s", not_connected_count)



if __name__ == "__main__":
    description = """
Determine connectivity for Hamming graphs based on taboo sets.

recommended usage:
cat input.txt | python run.py
"""

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument(
        '--organism-regex', default='^OS\s+(?P<name>[^\r]*)',
        )
    parser.add_argument(
        '--enzyme-regex', default='^ET\s+RM?2',
        )
    parser.add_argument(
        '--rsite-regex', default='^RS\s+(?P<nucleotides>[^,]*),',
        )
    parser.add_argument(
        '--max-wildcard', default=float('inf'), type=int,
        )
    parser.add_argument(
        '--log-level', default='DEBUG',
        )
    args = parser.parse_args()

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    ch.setLevel(args.log_level)
    logger.addHandler(ch)

    run(args.organism_regex, args.enzyme_regex, args.rsite_regex, args.max_wildcard)
