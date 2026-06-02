# Runtime Monitoring

This pipeline is intended to run with Docker on the study project VM. Use the Docker profile for normal runs:

```bash
nextflow run Anton-Bch/isoform-nf-core \
    -profile docker \
    --input samplesheet.csv \
    --outdir results
```

For test runs, the repository defaults to `test,docker` in `nf-test.config`, so `nf-test test` exercises the containerised runtime instead of expecting tools such as FastQC to exist on the host.

## VM-local monitor

Use `bin/nextflow-vm-monitor` to inspect active and recent Nextflow runs across a VM:

```bash
bin/nextflow-vm-monitor --root /home/ubuntu
```

It combines two sources:

- live Nextflow driver processes discovered from Linux `/proc`
- recent runs discovered from `.nextflow/history` files below the scan root

Useful variants:

```bash
# Refresh every 30 seconds
bin/nextflow-vm-monitor --root /home/ubuntu --watch 30

# JSON output for dashboards or cron jobs
bin/nextflow-vm-monitor --root /home/ubuntu --json

# Scan more deeply if projects are nested
bin/nextflow-vm-monitor --root /home/ubuntu --max-depth 6
```

This monitor does not require a service, database, or network connection. It is meant as the lightweight baseline for "what is running on this VM right now?" Nextflow also writes per-run reports, timelines, traces, DAGs, and execution history; see the [Nextflow reports documentation](https://docs.seqera.io/nextflow/reports) for the built-in reporting layer.

## Existing monitoring options

For production-grade monitoring, [Seqera Platform](https://docs.seqera.io/platform-cloud/getting-started/deployment-options) is the official Nextflow/Seqera solution. It gives a web UI for runs, task status, logs, resource usage, and team collaboration, but requires connecting the VM runs to a Seqera workspace.

For custom local infrastructure, Nextflow also has the [`nf-weblog`](https://registry.nextflow.io/plugins/nf-weblog) plugin. It can send workflow and task lifecycle events to an HTTP endpoint. A good next step would be a tiny local receiver service that stores those events in SQLite and serves a VM dashboard. That gives richer real-time state than scanning `.nextflow/history`, while still keeping the deployment local to the VM.
