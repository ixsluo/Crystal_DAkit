import click


@click.group()
@click.pass_context
def cli(ctx):
    pass


@cli.command()
def main():
    pass


if __name__ == "__main__":
    cli()

