from unidesign import ComputeStabilityConfig, UniDesignRunner, discover_binary
from unidesign.jobs import StabilityComputationJob
from unidesign.paths import project_root
import re


runner = UniDesignRunner(discover_binary())
pdb_path = project_root() / "example/ProteinProteinInteractionDesign/1ay7/1ay7.pdb"

job = StabilityComputationJob(runner, ComputeStabilityConfig(pdb_path=pdb_path))
result = job.run()
# Extract stdout text
stdout_text = result.run.stdout
#print("Raw output:\n", stdout_text)

# Parse the "Total" energy line
match = re.search(r"Total\s*=\s*([-0-9.]+)", stdout_text)
if match:
    total_energy = float(match.group(1))
    print("Total energy:", total_energy)

result.close()

