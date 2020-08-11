cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "medcelerate/farnsworth:latest"
  InlineJavascriptRequirement: {}

inputs:
  InputFiles:
    type: File[]
    inputBinding:
      position: 100

  isGenReion:
    type: boolean?
    inputBinding:
      prefix: "--gen_region"

  Output:
    type: string
    default: "consensus.vcf"
    inputBinding:
      prefix: "--output"
      valueFrom: "consensus.vcf"

baseCommand: ["/bin/Farnsworth"]

outputs:
  consensus_vcf:
    type: File
    outputBinding:
      glob: $(inputs.Output)

  regions:
    type: ["null", File]
    outputBinding:
      glob: "regions.txt"
