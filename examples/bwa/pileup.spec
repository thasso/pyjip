{
    "pileup": {
        "threads": 1, 
        "queue": "himem", 
        "priority": null, 
        "time": 360, 
        "memory": 0, 
        "account": "FB", 
        "extra": [], 
        "jobs": {
            "bwa_sam": {
                "threads": 2
            }, 
            "bam_index": {
                "threads": 3
            }, 
            "duplicates": {}, 
            "sam2bam": {"threads": 7}, 
            "mpileup": {
                "env":{
                    "PATH": "/apps/SAMTOOLS/0.1.19/bin:${PATH}"
                }
            }, 
            "bwa_align": {}
        }
    }
}
