{
    "queue": "Test",
    "jobs":{
        "first": {
            "queue": "second_queue"
        }, 
        "other":{
            "time": "5h",
            "threads": 10,
            "env":{
                "PATH": "/more:${PATH}"
            }
        }
    }
}
