import multiprocessing
import sys
from threading import Thread
import click
import data_io
import pairing


# Todo Find out if reads are being dropped
# Todo support for un-paired reads


def process_template(read_template):
    data_io.sam_to_array(read_template)

    # if read_template["name"] == "simulated_reads.0.10-id217_A_chr21:46697360_B_chr6:157282148-28985":
    #     click.echo(read_template["read1_length"], err=True)
    #     click.echo(read_template["read2_length"], err=True)
    #     for item in read_template["data"].tolist():
    #         click.echo(item, err=True)

    res = pairing.process(read_template)
    if res:
        data_io.apply_filter(read_template, *res)


def worker(queue, out_queue):
    while True:
        job = queue.get(True)

        if job == "Done":
            queue.task_done()
            out_queue.put("Done")
            break

        else:
            big_string = ""
            for data_tuple in job:
                read_template = data_io.make_template(*data_tuple)
                process_template(read_template)
                if read_template['passed']:
                    outstring = data_io.to_output(read_template)
                    if outstring:
                        big_string += outstring
                    else:
                        pass  # No mappings

            if len(big_string) > 0:
                out_queue.put(big_string)
            else:
                click.echo("WARNING: no output from job.", err=True)
        queue.task_done()


def process_reads(args):

    if not args["include"]:
        args["bias"] = 1.0
    else:
        click.echo("Elevating alignments in --include with --bias {}".format(args["bias"]), err=True)

    if args["outsam"] == "-" or args["output"] is None:
        click.echo("Writing alignments to stdout", err=True)
        outsam = sys.stdout
    else:
        click.echo("Writing alignments to {}".format(args["outsam"]), err=True)
        outsam = open(args["outsam"], "w")

    count = 0

    # Use multiprocessing:
    # https://stackoverflow.com/questions/17241663/filling-a-queue-and-managing-multiprocessing-in-python
    if args["procs"] != 1:

        cpus = args["procs"] if args["procs"] != 0 else multiprocessing.cpu_count()
        click.echo("fufi align runnning {} cpus".format(cpus), err=True)

        the_queue = multiprocessing.JoinableQueue(maxsize=cpus+2)
        out_queue = multiprocessing.Queue()

        the_pool = multiprocessing.Pool(args["procs"] if args["procs"] != 0 else multiprocessing.cpu_count(),
                                        worker, (the_queue, out_queue,))

        def writer_thread(q, outsam):
            while True:
                aln = q.get()
                if aln == "Done":
                    break
                elif aln == "Job failed":
                    click.echo("job failed", err=True)
                elif len(aln) > 1:
                    outsam.write(aln)

            click.echo("Writing done", err=True)

        writer = Thread(target=writer_thread, args=(out_queue, outsam, ))
        writer.setDaemon(True)
        writer.start()

        job = []
        itr = data_io.iterate_mappings(args)
        header_string = next(itr)
        out_queue.put(header_string)

        for data_tuple in itr:
            count += 1

            job.append(data_tuple)
            if len(job) > 500:
                the_queue.put(job)
                job = []
        if len(job) > 0:
            the_queue.put(job)

        the_queue.join()  # Wait for jobs to finish
        the_queue.put("Done")  # Send message to stop workers
        writer.join()  # Wait for writer to closing

    # Use single process for debugging
    else:
        click.echo("Single process", err=True)
        to_write = []  # Batch write

        itr = data_io.iterate_mappings(args)
        header_string = next(itr)
        to_write.append(header_string)

        for data_tuple in itr:

            count += 1
            temp = data_io.make_template(*data_tuple)
            process_template(temp)

            if temp['passed']:
                to_write.append(data_io.to_output(temp))

            if len(to_write) > 10000:  # Alignments to write
                for item in to_write:

                    if item is not None:
                        outsam.write(item)

                to_write = []

        for item in to_write:
            if item:
                outsam.write(item)

    if args["outsam"] != "-" or args["output"] is not None:
        outsam.close()
